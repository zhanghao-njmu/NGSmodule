/**
 * Task List Page - Task monitoring with real-time updates
 */
import React, { useEffect, useState } from 'react'
import {
  Button,
  Tag,
  Progress,
  Space,
  Select,
} from 'antd'
import {
  PlusOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ClockCircleOutlined,
  StopOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useTaskStore } from '../../store/taskStore'
import { useProjectStore } from '../../store/projectStore'
import { websocketService } from '../../services/websocket.service'
import { PageHeader, DataTable, StatisticCard, StatusTag } from '../../components/common'
import { confirm } from '../../components/common/ConfirmDialog'
import { toast } from '../../utils/notification'
import type { StatisticItem } from '../../components/common'
import type { Task, TaskStatus } from '../../types/task'
import dayjs from 'dayjs'

const { Option } = Select

export const TaskList: React.FC = () => {
  const {
    tasks,
    stats,
    loading,
    fetchTasks,
    fetchStats,
    cancelTask,
    handleWebSocketMessage,
  } = useTaskStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')

  useEffect(() => {
    fetchProjects()
    fetchStats()

    // Setup WebSocket for real-time updates
    const token = localStorage.getItem('auth_token')
    if (token) {
      websocketService.connect(token)
      websocketService.addMessageHandler(handleWebSocketMessage)
    }

    return () => {
      websocketService.disconnect()
    }
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchTasks({ project_id: selectedProject })
    } else {
      fetchTasks()
    }
  }, [selectedProject])

  const handleCancelTask = (task: Task) => {
    confirm({
      title: '取消任务',
      content: `您确定要取消任务 "${task.task_name}" 吗？`,
      okText: '取消任务',
      cancelText: '保持运行',
      okButtonProps: { danger: true },
      onConfirm: async () => {
        const loadingToast = toast.loading('取消中...')
        try {
          await cancelTask(task.id)
          loadingToast()
          toast.success('任务已取消')
          fetchTasks()
          fetchStats()
        } catch (error) {
          loadingToast()
          toast.error('取消失败，请重试')
        }
      }
    })
  }

  const columns: ColumnsType<Task> = [
    {
      title: 'Task Name',
      dataIndex: 'task_name',
      key: 'task_name',
      render: (name) => <span style={{ fontWeight: 500 }}>{name}</span>,
    },
    {
      title: 'Type',
      dataIndex: 'task_type',
      key: 'task_type',
      width: 120,
      render: (type) => type && <Tag color="blue">{type}</Tag>,
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 120,
      render: (status: TaskStatus) => <StatusTag status={status} />,
    },
    {
      title: 'Progress',
      dataIndex: 'progress',
      key: 'progress',
      width: 150,
      render: (progress) => (
        <Progress
          percent={Math.round(progress)}
          size="small"
          status={progress === 100 ? 'success' : 'active'}
        />
      ),
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 100,
      render: (_, record) => (
        <Space>
          {record.status === 'running' && (
            <Button
              size="small"
              danger
              icon={<StopOutlined />}
              onClick={() => handleCancelTask(record)}
            >
              Cancel
            </Button>
          )}
        </Space>
      ),
    },
  ]

  const statisticItems: StatisticItem[] = [
    {
      key: 'total',
      title: 'Total',
      value: stats?.total_tasks || 0,
      valueStyle: { color: 'var(--color-primary)' },
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'running',
      title: 'Running',
      value: stats?.running_tasks || 0,
      valueStyle: { color: 'var(--color-warning)' },
      prefix: <SyncOutlined spin />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'completed',
      title: 'Completed',
      value: stats?.completed_tasks || 0,
      valueStyle: { color: 'var(--color-success)' },
      prefix: <CheckCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'failed',
      title: 'Failed',
      value: stats?.failed_tasks || 0,
      valueStyle: { color: 'var(--color-error)' },
      prefix: <CloseCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
  ]

  return (
    <div>
      {/* Statistics */}
      <StatisticCard items={statisticItems} />

      <PageHeader
        left={
          <Select
            placeholder="All Projects"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
            allowClear
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
        }
        right={
          <Button type="primary" icon={<PlusOutlined />}>
            New Task
          </Button>
        }
      />

      <DataTable
        columns={columns}
        dataSource={tasks}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Tasks"
        emptyDescription="No tasks have been created yet"
      />
    </div>
  )
}

export default TaskList
