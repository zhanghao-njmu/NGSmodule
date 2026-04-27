/**
 * Admin Dashboard - System overview and user management
 */
import type React from 'react'
import { useState } from 'react'
import {
  Card,
  Tag,
  Button,
  Space,
  Modal,
  Form,
  Input,
  Select,
  InputNumber,
  Switch,
  Progress,
  Popconfirm,
  Typography,
} from 'antd'
import {
  UserOutlined,
  ProjectOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  DatabaseOutlined,
  EditOutlined,
  DeleteOutlined,
  StopOutlined,
  PlayCircleOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'

import {
  useSystemStats,
  useAdminUsers,
  useUpdateAdminUser,
  useToggleUserStatus,
  useDeleteAdminUser,
} from '@/hooks/queries'
import {
  StatisticCard,
  DataTable,
  StatusTag,
  PageSkeleton,
  FadeIn,
  StaggeredList,
  EnhancedEmptyState,
} from '@/components/common'
import { toast } from '@/utils/notification'
import { formatFileSize } from '@/utils/format'
import type { StatisticItem } from '@/components/common'
import type { User, UserAdminUpdate, SystemStats } from '@/types/admin'
import dayjs from 'dayjs'

const { Option } = Select
const { Title, Text } = Typography

export const AdminDashboard: React.FC = () => {
  const [editModalOpen, setEditModalOpen] = useState(false)
  const [selectedUser, setSelectedUser] = useState<User | null>(null)
  const [form] = Form.useForm()

  const { data: stats, isLoading: statsLoading } = useSystemStats() as {
    data: SystemStats | undefined
    isLoading: boolean
  }
  const { data: usersData, isLoading: usersLoading } = useAdminUsers({ limit: 100 })
  const users: User[] = (usersData as any)?.users ?? (usersData as any) ?? []

  const updateMutation = useUpdateAdminUser()
  const toggleMutation = useToggleUserStatus()
  const deleteMutation = useDeleteAdminUser()

  const loading = statsLoading || usersLoading

  const handleEdit = (user: User) => {
    setSelectedUser(user)
    form.setFieldsValue({
      full_name: user.full_name,
      email: user.email,
      organization: user.organization,
      role: user.role,
      is_active: user.is_active,
      storage_quota: Math.round((user?.storage_quota ?? 0) / (1024 * 1024 * 1024)),
    })
    setEditModalOpen(true)
  }

  const handleEditSubmit = async () => {
    if (!selectedUser) {
      return
    }
    try {
      const values = await form.validateFields()
      const updateData: UserAdminUpdate = {
        ...values,
        storage_quota: values.storage_quota * 1024 * 1024 * 1024,
      }
      await updateMutation.mutateAsync({ id: selectedUser.id, data: updateData })
      toast.success('User updated successfully')
      setEditModalOpen(false)
    } catch (error: any) {
      toast.error(`Failed to update user: ${error.message}`)
    }
  }

  const handleToggleStatus = async (user: User) => {
    try {
      await toggleMutation.mutateAsync({ id: user.id, isActive: !user.is_active })
      toast.success(`User ${user.is_active ? 'deactivated' : 'activated'} successfully`)
    } catch (error: any) {
      toast.error(`Failed to toggle user status: ${error.message}`)
    }
  }

  const handleDelete = async (userId: string) => {
    try {
      await deleteMutation.mutateAsync({ id: userId })
      toast.success('User deleted successfully')
    } catch (error: any) {
      toast.error(`Failed to delete user: ${error.message}`)
    }
  }

  const columns: ColumnsType<User> = [
    {
      title: 'Username',
      dataIndex: 'username',
      key: 'username',
      render: (username) => <span style={{ fontWeight: 500 }}>{username}</span>,
    },
    {
      title: 'Email',
      dataIndex: 'email',
      key: 'email',
    },
    {
      title: 'Role',
      dataIndex: 'role',
      key: 'role',
      width: 100,
      render: (role) => <Tag color={role === 'admin' ? 'red' : 'blue'}>{role.toUpperCase()}</Tag>,
    },
    {
      title: 'Status',
      dataIndex: 'is_active',
      key: 'is_active',
      width: 100,
      render: (is_active) => <StatusTag status={is_active ? 'active' : 'inactive'} />,
    },
    {
      title: 'Storage',
      key: 'storage',
      width: 200,
      render: (_, record) => {
        const percent = ((record?.storage_used ?? 0) / (record?.storage_quota ?? 1)) * 100
        return (
          <div>
            <Progress percent={Math.round(percent)} size="small" status={percent > 90 ? 'exception' : 'normal'} />
            <div style={{ fontSize: 12, color: '#666' }}>
              {formatFileSize(record?.storage_used ?? 0)} / {formatFileSize(record?.storage_quota ?? 0)}
            </div>
          </div>
        )
      },
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 120,
      render: (date) => dayjs(date).format('YYYY-MM-DD'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 180,
      render: (_, record) => (
        <Space>
          <Button size="small" icon={<EditOutlined />} onClick={() => handleEdit(record)}>
            Edit
          </Button>
          <Button
            size="small"
            icon={record.is_active ? <StopOutlined /> : <PlayCircleOutlined />}
            onClick={() => handleToggleStatus(record)}
          >
            {record.is_active ? 'Deactivate' : 'Activate'}
          </Button>
          <Popconfirm
            title="Delete user"
            description="Are you sure? This will delete all user data."
            onConfirm={() => handleDelete(record.id)}
            okText="Yes"
            cancelText="No"
          >
            <Button size="small" danger icon={<DeleteOutlined />} />
          </Popconfirm>
        </Space>
      ),
    },
  ]

  const statisticItems: StatisticItem[] = [
    {
      key: 'users',
      title: 'Total Users',
      value: stats?.total_users || 0,
      prefix: <UserOutlined />,
      valueStyle: { color: 'var(--color-primary)' },
      suffix: (
        <div style={{ fontSize: 12, color: 'var(--text-secondary)', marginTop: 8 }}>
          {stats?.active_users || 0} active
        </div>
      ),
    },
    {
      key: 'items',
      title: 'Total Projects',
      value: stats?.total_projects || 0,
      prefix: <ProjectOutlined />,
      valueStyle: { color: 'var(--color-success)' },
    },
    {
      key: 'samples',
      title: 'Total Samples',
      value: stats?.total_samples || 0,
      prefix: <ExperimentOutlined />,
      valueStyle: { color: 'var(--color-warning)' },
    },
    {
      key: 'tasks',
      title: 'Total Tasks',
      value: stats?.total_tasks || 0,
      prefix: <ThunderboltOutlined />,
      valueStyle: { color: '#722ed1' },
      suffix: (
        <div style={{ fontSize: 12, marginTop: 8 }}>
          <CheckCircleOutlined style={{ color: 'var(--color-success)' }} /> {stats?.completed_tasks || 0} |{' '}
          <CloseCircleOutlined style={{ color: 'var(--color-error)' }} /> {stats?.failed_tasks || 0}
        </div>
      ),
    },
  ]

  if (loading && !stats && users.length === 0) {
    return <PageSkeleton hasHeader rows={8} />
  }

  return (
    <div>
      {/* Header with animation */}
      <FadeIn direction="up" delay={0} duration={300}>
        <Space align="center" style={{ marginBottom: 24 }}>
          <UserOutlined style={{ fontSize: 28, color: 'var(--color-primary)' }} />
          <div>
            <Title level={2} style={{ margin: 0 }}>
              Admin Dashboard
            </Title>
            <Text type="secondary">System overview and user management</Text>
          </div>
        </Space>
      </FadeIn>

      {/* System Statistics with staggered animation */}
      <StaggeredList staggerDelay={80} baseDelay={0} direction="up">
        <StatisticCard items={statisticItems} gutter={[16, 16]} />
      </StaggeredList>

      {/* Storage Statistics with animation */}
      <FadeIn direction="up" delay={100} duration={300}>
        <Card
          title={
            <>
              <DatabaseOutlined /> Storage Usage
            </>
          }
          style={{ marginTop: 24, marginBottom: 24 }}
        >
          <Progress
            percent={Math.round(((stats?.total_storage_used || 0) / (stats?.total_storage_quota || 1)) * 100)}
            status="active"
          />
          <div style={{ marginTop: 8, fontSize: 14 }}>
            {formatFileSize(stats?.total_storage_used || 0)} / {formatFileSize(stats?.total_storage_quota || 0)}
          </div>
        </Card>
      </FadeIn>

      {/* User Management Table with animation */}
      <FadeIn direction="up" delay={200} duration={300}>
        <DataTable
          title="User Management"
          columns={columns}
          dataSource={users}
          rowKey="id"
          loading={loading}
          pagination={{ pageSize: 20, showSizeChanger: true, showTotal: (total) => `Total ${total} users` }}
          locale={{
            emptyText: (
              <EnhancedEmptyState
                type="noData"
                title="No users yet"
                description="No users have been registered yet. Users will appear here after registration."
                size="default"
              />
            ),
          }}
        />
      </FadeIn>

      {/* Edit User Modal */}
      <Modal
        title="Edit User"
        open={editModalOpen}
        onCancel={() => setEditModalOpen(false)}
        onOk={handleEditSubmit}
        width={600}
        confirmLoading={updateMutation.isPending}
      >
        <Form form={form} layout="vertical">
          <Form.Item name="full_name" label="Full Name">
            <Input placeholder="Enter full name" />
          </Form.Item>

          <Form.Item name="email" label="Email" rules={[{ type: 'email' }]}>
            <Input placeholder="Enter email" />
          </Form.Item>

          <Form.Item name="organization" label="Organization">
            <Input placeholder="Enter organization" />
          </Form.Item>

          <Form.Item name="role" label="Role" rules={[{ required: true }]}>
            <Select>
              <Option value="user">User</Option>
              <Option value="admin">Admin</Option>
            </Select>
          </Form.Item>

          <Form.Item name="is_active" label="Active" valuePropName="checked">
            <Switch />
          </Form.Item>

          <Form.Item
            name="storage_quota"
            label="Storage Quota (GB)"
            rules={[{ required: true, type: 'number', min: 1 }]}
          >
            <InputNumber style={{ width: '100%' }} min={1} />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  )
}

export default AdminDashboard
