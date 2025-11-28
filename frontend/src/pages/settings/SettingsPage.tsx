import { useState, useEffect, useCallback } from 'react'
import {
  Card,
  Row,
  Col,
  Typography,
  Form,
  Input,
  Switch,
  Button,
  Select,
  Divider,
  Space,
  message,
  Modal,
  Table,
  Tag,
  Tooltip,
  Alert,
} from 'antd'
import {
  SettingOutlined,
  BellOutlined,
  LockOutlined,
  GlobalOutlined,
  ApiOutlined,
  DeleteOutlined,
  PlusOutlined,
  CopyOutlined,
  EyeOutlined,
  EyeInvisibleOutlined,
  SafetyOutlined,
  ThunderboltOutlined,
  MailOutlined,
} from '@ant-design/icons'
import { PageHeader } from '@/components/common'
import { useTheme } from '@/store/themeStore'
import { userService, type ApiToken } from '@/services/user.service'
import { DesignTokens } from '@/styles/design-tokens'
import { logger } from '@/utils/logger'
import './SettingsPage.css'

const { Text, Paragraph } = Typography
const { Option } = Select

export const SettingsPage: React.FC = () => {
  const { mode, toggleMode } = useTheme()
  const [accountForm] = Form.useForm()
  const [notificationForm] = Form.useForm()
  const [privacyForm] = Form.useForm()
  const [newTokenModalVisible, setNewTokenModalVisible] = useState(false)
  const [tokenForm] = Form.useForm()
  const [visibleTokens, setVisibleTokens] = useState<Set<string>>(new Set())
  const [apiTokens, setApiTokens] = useState<ApiToken[]>([])
  const [loading, setLoading] = useState(false)

  // Fetch API tokens on mount
  const fetchApiTokens = useCallback(async () => {
    try {
      const tokens = await userService.getApiTokens()
      setApiTokens(tokens)
    } catch (error) {
      logger.error('Failed to fetch API tokens:', error)
    }
  }, [])

  useEffect(() => {
    fetchApiTokens()
  }, [fetchApiTokens])

  // Account Settings
  const handleAccountSettingsSave = async (values: Record<string, unknown>) => {
    setLoading(true)
    try {
      await userService.updateSettings({
        language: values.language as string,
        timezone: values.timezone as string,
        dateFormat: values.dateFormat as string,
        theme: mode,
      })
      message.success('Account settings saved successfully!')
    } catch (error) {
      logger.error('Failed to save account settings:', error)
      message.error('Failed to save account settings')
    } finally {
      setLoading(false)
    }
  }

  // Notification Settings
  const handleNotificationSettingsSave = async (values: Record<string, unknown>) => {
    setLoading(true)
    try {
      await userService.updateNotificationSettings({
        emailNotifications: values.emailNotifications as boolean,
        pipelineComplete: values.pipelineComplete as boolean,
        taskFailed: values.taskFailed as boolean,
        systemUpdates: values.systemUpdates as boolean,
        weeklyReport: values.weeklyReport as boolean,
        browserNotifications: values.browserNotifications as boolean,
      })
      message.success('Notification preferences saved!')
    } catch (error) {
      logger.error('Failed to save notification settings:', error)
      message.error('Failed to save notification settings')
    } finally {
      setLoading(false)
    }
  }

  // Privacy Settings
  const handlePrivacySettingsSave = async (values: Record<string, unknown>) => {
    setLoading(true)
    try {
      await userService.updatePrivacySettings({
        profileVisible: values.profileVisible as boolean,
        activityVisible: values.activityVisible as boolean,
        twoFactorAuth: values.twoFactorAuth as boolean,
        sessionTimeout: values.sessionTimeout as number,
      })
      message.success('Privacy settings updated!')
    } catch (error) {
      logger.error('Failed to update privacy settings:', error)
      message.error('Failed to update privacy settings')
    } finally {
      setLoading(false)
    }
  }

  // API Token Management
  const handleCreateToken = async (values: { name: string; description?: string }) => {
    try {
      const newToken = await userService.createApiToken({
        name: values.name,
        description: values.description,
      })
      setApiTokens([...apiTokens, newToken])
      message.success('API token created successfully!')
      setNewTokenModalVisible(false)
      tokenForm.resetFields()

      // Show full token once
      Modal.info({
        title: 'New API Token Created',
        content: (
          <div>
            <Alert
              message="Save this token now!"
              description="For security reasons, you won't be able to see this token again. Please copy and save it securely."
              type="warning"
              showIcon
              style={{ marginBottom: 16 }}
            />
            <Paragraph copyable={{ text: newToken.token }}>
              <Text code style={{ fontSize: 12 }}>
                {newToken.token}
              </Text>
            </Paragraph>
          </div>
        ),
        okText: 'I have saved the token',
      })
    } catch (error) {
      logger.error('Failed to create API token:', error)
      message.error('Failed to create API token')
    }
  }

  const handleDeleteToken = (tokenId: string) => {
    Modal.confirm({
      title: 'Delete API Token',
      content: 'Are you sure you want to delete this token? This action cannot be undone.',
      okText: 'Delete',
      okType: 'danger',
      cancelText: 'Cancel',
      onOk: async () => {
        try {
          await userService.deleteApiToken(tokenId)
          setApiTokens(apiTokens.filter((token) => token.id !== tokenId))
          message.success('API token deleted successfully')
        } catch (error) {
          logger.error('Failed to delete API token:', error)
          message.error('Failed to delete API token')
        }
      },
    })
  }

  const toggleTokenVisibility = (tokenId: string) => {
    const newVisible = new Set(visibleTokens)
    if (newVisible.has(tokenId)) {
      newVisible.delete(tokenId)
    } else {
      newVisible.add(tokenId)
    }
    setVisibleTokens(newVisible)
  }

  const handleCopyToken = (token: string) => {
    navigator.clipboard.writeText(token)
    message.success('Token copied to clipboard')
  }

  const tokenColumns = [
    {
      title: 'Name',
      dataIndex: 'name',
      key: 'name',
    },
    {
      title: 'Token',
      dataIndex: 'token',
      key: 'token',
      render: (token: string, record: ApiToken) => {
        const isVisible = visibleTokens.has(record.id)
        return (
          <Space>
            <Text code style={{ fontSize: 12 }}>
              {isVisible ? token : '••••••••••••••••'}
            </Text>
            <Tooltip title={isVisible ? 'Hide' : 'Show'}>
              <Button
                type="text"
                size="small"
                icon={isVisible ? <EyeInvisibleOutlined /> : <EyeOutlined />}
                onClick={() => toggleTokenVisibility(record.id)}
              />
            </Tooltip>
            <Tooltip title="Copy">
              <Button type="text" size="small" icon={<CopyOutlined />} onClick={() => handleCopyToken(token)} />
            </Tooltip>
          </Space>
        )
      },
    },
    {
      title: 'Created',
      dataIndex: 'createdAt',
      key: 'createdAt',
    },
    {
      title: 'Last Used',
      dataIndex: 'lastUsed',
      key: 'lastUsed',
      render: (lastUsed?: string) => lastUsed || 'Never',
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      render: (status: string) => <Tag color={status === 'active' ? 'success' : 'error'}>{status.toUpperCase()}</Tag>,
    },
    {
      title: 'Action',
      key: 'action',
      render: (_: any, record: ApiToken) => (
        <Button type="text" danger size="small" icon={<DeleteOutlined />} onClick={() => handleDeleteToken(record.id)}>
          Delete
        </Button>
      ),
    },
  ]

  return (
    <div className="settings-page">
      <PageHeader title="设置" subtitle="管理您的偏好和应用程序配置" icon={<SettingOutlined />} />

      <Row gutter={[24, 24]}>
        {/* Account Settings */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <GlobalOutlined style={{ color: DesignTokens.colors.primary.main }} />
                <Text strong>账户设置</Text>
              </Space>
            }
            className="settings-card"
          >
            <Form
              form={accountForm}
              layout="vertical"
              onFinish={handleAccountSettingsSave}
              initialValues={{
                language: 'zh-CN',
                timezone: 'Asia/Shanghai',
                dateFormat: 'YYYY-MM-DD',
                theme: mode,
              }}
            >
              <Form.Item label="语言" name="language">
                <Select>
                  <Option value="zh-CN">简体中文</Option>
                  <Option value="en-US">English</Option>
                  <Option value="ja-JP">日本語</Option>
                </Select>
              </Form.Item>

              <Form.Item label="时区" name="timezone">
                <Select showSearch>
                  <Option value="Asia/Shanghai">Asia/Shanghai (UTC+8)</Option>
                  <Option value="America/New_York">America/New_York (UTC-5)</Option>
                  <Option value="Europe/London">Europe/London (UTC+0)</Option>
                  <Option value="Asia/Tokyo">Asia/Tokyo (UTC+9)</Option>
                </Select>
              </Form.Item>

              <Form.Item label="日期格式" name="dateFormat">
                <Select>
                  <Option value="YYYY-MM-DD">YYYY-MM-DD</Option>
                  <Option value="MM/DD/YYYY">MM/DD/YYYY</Option>
                  <Option value="DD/MM/YYYY">DD/MM/YYYY</Option>
                </Select>
              </Form.Item>

              <Form.Item label="主题模式">
                <Space>
                  <Switch
                    checked={mode === 'dark'}
                    onChange={toggleMode}
                    checkedChildren="暗色"
                    unCheckedChildren="亮色"
                  />
                  <Text type="secondary">{mode === 'dark' ? '暗色模式已启用' : '亮色模式已启用'}</Text>
                </Space>
              </Form.Item>

              <Form.Item>
                <Button type="primary" htmlType="submit" block loading={loading}>
                  保存设置
                </Button>
              </Form.Item>
            </Form>
          </Card>
        </Col>

        {/* Notification Settings */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <BellOutlined style={{ color: DesignTokens.colors.warning.main }} />
                <Text strong>通知设置</Text>
              </Space>
            }
            className="settings-card"
          >
            <Form
              form={notificationForm}
              layout="vertical"
              onFinish={handleNotificationSettingsSave}
              initialValues={{
                emailNotifications: true,
                pipelineComplete: true,
                taskFailed: true,
                systemUpdates: false,
                weeklyReport: true,
                browserNotifications: true,
              }}
            >
              <Paragraph type="secondary" style={{ marginBottom: 16 }}>
                选择您希望接收的通知类型
              </Paragraph>

              <Form.Item name="emailNotifications" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>
                      <MailOutlined /> 邮件通知
                    </Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      通过邮件接收重要通知
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="pipelineComplete" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>
                      <ThunderboltOutlined /> 流程完成
                    </Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      当分析流程完成时通知我
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="taskFailed" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>任务失败</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      当任务失败时立即通知我
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="weeklyReport" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>每周报告</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      接收每周活动摘要
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="browserNotifications" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>浏览器通知</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      启用浏览器推送通知
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Form.Item>
                <Button type="primary" htmlType="submit" block loading={loading}>
                  保存通知设置
                </Button>
              </Form.Item>
            </Form>
          </Card>
        </Col>

        {/* Privacy & Security */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <SafetyOutlined style={{ color: DesignTokens.colors.success.main }} />
                <Text strong>隐私与安全</Text>
              </Space>
            }
            className="settings-card"
          >
            <Form
              form={privacyForm}
              layout="vertical"
              onFinish={handlePrivacySettingsSave}
              initialValues={{
                profileVisible: true,
                activityVisible: false,
                twoFactorAuth: false,
                sessionTimeout: 30,
              }}
            >
              <Form.Item name="profileVisible" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>公开个人资料</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      允许其他用户查看您的个人资料
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="activityVisible" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>显示活动状态</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      让其他人看到您的在线状态
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item name="twoFactorAuth" valuePropName="checked">
                <div className="switch-item">
                  <div>
                    <Text strong>
                      <LockOutlined /> 双因素认证
                    </Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      为您的账户添加额外的安全层
                    </Text>
                  </div>
                  <Switch />
                </div>
              </Form.Item>

              <Divider style={{ margin: '12px 0' }} />

              <Form.Item label="会话超时（分钟）" name="sessionTimeout">
                <Select>
                  <Option value={15}>15 分钟</Option>
                  <Option value={30}>30 分钟</Option>
                  <Option value={60}>1 小时</Option>
                  <Option value={120}>2 小时</Option>
                  <Option value={0}>永不超时</Option>
                </Select>
              </Form.Item>

              <Form.Item>
                <Button type="primary" htmlType="submit" block loading={loading}>
                  保存隐私设置
                </Button>
              </Form.Item>
            </Form>
          </Card>
        </Col>

        {/* API Token Management */}
        <Col xs={24}>
          <Card
            title={
              <Space>
                <ApiOutlined style={{ color: DesignTokens.colors.info.main }} />
                <Text strong>API 令牌管理</Text>
              </Space>
            }
            extra={
              <Button type="primary" icon={<PlusOutlined />} onClick={() => setNewTokenModalVisible(true)}>
                创建新令牌
              </Button>
            }
            className="settings-card"
          >
            <Alert
              message="API 令牌安全提示"
              description="API 令牌允许您通过编程方式访问 NGSmodule。请妥善保管您的令牌，不要与他人分享。"
              type="info"
              showIcon
              style={{ marginBottom: 16 }}
            />

            <Table columns={tokenColumns} dataSource={apiTokens} rowKey="id" pagination={false} />
          </Card>
        </Col>
      </Row>

      {/* Create Token Modal */}
      <Modal
        title="创建新的 API 令牌"
        open={newTokenModalVisible}
        onCancel={() => {
          setNewTokenModalVisible(false)
          tokenForm.resetFields()
        }}
        onOk={() => tokenForm.submit()}
        okText="创建"
        cancelText="取消"
      >
        <Form form={tokenForm} layout="vertical" onFinish={handleCreateToken}>
          <Form.Item label="令牌名称" name="name" rules={[{ required: true, message: '请输入令牌名称' }]}>
            <Input placeholder="例如：Pipeline Automation" />
          </Form.Item>

          <Form.Item label="描述（可选）" name="description">
            <Input.TextArea rows={3} placeholder="描述此令牌的用途..." />
          </Form.Item>

          <Alert message="创建后您将只能看到一次完整的令牌，请务必保存" type="warning" showIcon />
        </Form>
      </Modal>
    </div>
  )
}
